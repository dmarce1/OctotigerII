/** @file
 * @brief Heterogeneous quantity tuples and componentwise state operations.
 * @ingroup units
 */
#pragma once
#include <array>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <utility>
#include "octotigerII/buildConfig.hpp"
#include "octotigerII/units/cgs.hpp"


namespace octotigerII::units {


// A state is a heterogeneous tuple of quantities. There is no raw numeric
// indexing. Runtime indexing is available only when every unit is identical.
/// Tuple of quantities with independently checked physical dimensions.
/// get<I>() is compile-time indexing; operator[] exists only for homogeneous types.
/// Named accessors follow the field ordering declared by each physics adapter.
/// Arithmetic acts componentwise; multiplying by a quantity changes all dimensions.
/// @ingroup units
template <typename... Q>
class State {
	std::tuple<Q...> data_{};

	template <typename F, std::size_t... I>
	constexpr void each(F&& f, std::index_sequence<I...>) {
		(f(std::integral_constant<std::size_t, I>{}, std::get<I>(data_)), ...);
	}

	template <typename F, std::size_t... I>
	constexpr void each(F&& f, std::index_sequence<I...>) const {
		(f(std::integral_constant<std::size_t, I>{}, std::get<I>(data_)), ...);
	}
	using First = std::tuple_element_t<0, std::tuple<Q...>>;

public:

	constexpr State() = default;

	constexpr State(Q... q)
	  : data_(q...) {}

	constexpr First& operator[](int i)
		requires((std::is_same_v<First, Q>) && ...)
	{
		auto pointers = std::apply([](auto&... q) { return std::array<First*, sizeof...(Q)>{&q...}; }, data_);
		return *pointers.at(i);
	}

	constexpr First const& operator[](int i) const
		requires((std::is_same_v<First, Q>) && ...)
	{
		auto pointers = std::apply([](auto const&... q) { return std::array<First const*, sizeof...(Q)>{&q...}; }, data_);
		return *pointers.at(i);
	}

	constexpr State& operator+=(State const& b) {
		forEach([&](auto i, auto& q) { q += b.template get<i>(); });
		return *this;
	}

	constexpr State& operator-=(State const& b) {
		forEach([&](auto i, auto& q) { q -= b.template get<i>(); });
		return *this;
	}

	constexpr State& operator*=(Real s) {
		forEach([&](auto, auto& q) { q *= s; });
		return *this;
	}

	constexpr State& operator/=(Real s) {
		forEach([&](auto, auto& q) { q /= s; });
		return *this;
	}

	/// Return the component selected by its compile-time tuple index.
	template <std::size_t I>
	constexpr auto& get() {
		return std::get<I>(data_);
	}

	/// Return the component selected by its compile-time tuple index.
	template <std::size_t I>
	constexpr auto const& get() const {
		return std::get<I>(data_);
	}

	/// Visit each component with its compile-time index and typed value.
	template <typename F>
	constexpr void forEach(F&& f) {
		each(std::forward<F>(f), std::index_sequence_for<Q...>{});
	}

	/// Visit each component with its compile-time index and typed value.
	template <typename F>
	constexpr void forEach(F&& f) const {
		each(std::forward<F>(f), std::index_sequence_for<Q...>{});
	}

	// Named views for Euler states and fluxes; each returns the actual stored
	// quantity, never a proxy or a number with units attached afterwards.
	/// Access component 0 of a gravity state.
	constexpr auto& potential() {
		return get<0>();
	}

	/// Access component 0 of a gravity state.
	constexpr auto const& potential() const {
		return get<0>();
	}

	/// Access component 1+axis of a gravity state.
	constexpr auto& acceleration(int axis) {
		return momentum(axis);
	}

	/// Access component 1+axis of a gravity state.
	constexpr auto const& acceleration(int axis) const {
		return momentum(axis);
	}

	/// Access the first component, interpreted as mass density by the hydro schema.
	constexpr auto& density() {
		return get<0>();
	}

	/// Access the first component, interpreted as mass density by the hydro schema.
	constexpr auto const& density() const {
		return get<0>();
	}

	/// Access the first component of a hydro flux (mass flux).
	constexpr auto& mass() {
		return get<0>();
	}

	/// Access the first component of a hydro flux (mass flux).
	constexpr auto const& mass() const {
		return get<0>();
	}

	/// Access component 1+axis; the hydro schema uses Cartesian momentum density.
	constexpr auto& momentum(int axis) {
		if (axis < 0 || axis >= ndim) throw std::out_of_range("State vector axis outside ndim");
		return [&]<std::size_t... I>(std::index_sequence<I...>) -> decltype(auto) {
			auto pointers = std::array{&get<I + 1>()...};
			return *pointers[axis];
		}(std::make_index_sequence<ndim>{});
	}

	/// Access component 1+axis; the hydro schema uses Cartesian momentum density.
	constexpr auto const& momentum(int axis) const {
		if (axis < 0 || axis >= ndim) throw std::out_of_range("State vector axis outside ndim");
		return [&]<std::size_t... I>(std::index_sequence<I...>) -> decltype(auto) {
			auto pointers = std::array{&get<I + 1>()...};
			return *pointers[axis];
		}(std::make_index_sequence<ndim>{});
	}

	/// Access component 1+axis of a primitive hydro state.
	constexpr auto& velocity(int axis) {
		return momentum(axis);
	}

	/// Access component 1+axis of a primitive hydro state.
	constexpr auto const& velocity(int axis) const {
		return momentum(axis);
	}

	constexpr auto& totalEnergy() {
		return get<ndim + 1>();
	}

	constexpr auto const& totalEnergy() const {
		return get<ndim + 1>();
	}

	/// Access component ndim+1 of a primitive hydro state.
	constexpr auto& pressure() {
		return get<ndim + 1>();
	}

	/// Access component ndim+1 of a primitive hydro state.
	constexpr auto const& pressure() const {
		return get<ndim + 1>();
	}

	/// Access component ndim+1 for a hydro state, otherwise component 0.
	constexpr auto& energy() {
		if constexpr (sizeof...(Q) == ndim + 2)
			return get<ndim + 1>();
		else
			return get<0>();
	}

	/// Access component ndim+1 for a hydro state, otherwise component 0.
	constexpr auto const& energy() const {
		if constexpr (sizeof...(Q) == ndim + 2)
			return get<ndim + 1>();
		else
			return get<0>();
	}

	/// Access component 1+axis; physical radiation states return EnergyFlux.
	constexpr auto& radiativeFlux(int axis) {
		return momentum(axis);
	}

	/// Access component 1+axis; physical radiation states return EnergyFlux.
	constexpr auto const& radiativeFlux(int axis) const {
		return momentum(axis);
	}

	constexpr void setDensity(First q) {
		get<0>() = q;
	}

	constexpr void setMass(First q) {
		get<0>() = q;
	}

	constexpr void setMomentum(int i, std::tuple_element_t<1, std::tuple<Q...>> q) {
		momentum(i) = q;
	}

	constexpr void setVelocity(int i, std::tuple_element_t<1, std::tuple<Q...>> q) {
		momentum(i) = q;
	}

	constexpr void setRadiativeFlux(int i, std::tuple_element_t<1, std::tuple<Q...>> q) {
		momentum(i) = q;
	}

	template <typename V>
	constexpr void setTotalEnergy(V q) {
		get<ndim + 1>() = q;
	}

	template <typename V>
	constexpr void setPressure(V q) {
		get<ndim + 1>() = q;
	}

	template <typename V>
	constexpr void setEnergy(V q) {
		energy() = q;
	}

	/// Serialize this value with its compile-time quantity types preserved.
	template <typename Archive>
	void serialize(Archive& a, unsigned) {
		forEach([&](auto, auto& q) { a & q; });
	}

	/// Return the number of represented elements or blocks.
	static constexpr int size() {
		return sizeof...(Q);
	}

	friend State componentAbs(State a) {
		a.forEach([](auto, auto& q) { q = units::abs(q); });
		return a;
	}

	friend bool finite(State const& a) {
		bool ok = true;
		a.forEach([&](auto, auto const& q) { ok = ok && units::finite(q); });
		return ok;
	}

	template <typename U>
	friend constexpr auto operator*(State const& a, boost::units::quantity<U, Real> s) {
		return std::apply([&](auto const&... q) { return State<decltype(q * s)...>(q * s...); }, a.data_);
	}

	template <typename U>
	friend constexpr auto operator*(boost::units::quantity<U, Real> s, State const& a) {
		return a * s;
	}

	template <typename U>
	friend constexpr auto operator/(State const& a, boost::units::quantity<U, Real> s) {
		return std::apply([&](auto const&... q) { return State<decltype(q / s)...>(q / s...); }, a.data_);
	}

	friend constexpr State operator+(State a, State const& b) {
		return a += b;
	}

	friend constexpr State operator-(State a, State const& b) {
		return a -= b;
	}

	friend constexpr State operator-(State a) {
		return a *= Real(-1);
	}

	friend constexpr State operator*(State a, Real s) {
		return a *= s;
	}

	friend constexpr State operator*(Real s, State a) {
		return a *= s;
	}

	friend constexpr State operator/(State a, Real s) {
		return a /= s;
	}
};


namespace detail {
template <typename Scalar, typename Vector, std::size_t... I>
auto scalarVector(std::index_sequence<I...>) -> State<Scalar, decltype((void) I, Vector{})...>;

template <typename Density, typename Vector, typename Energy, std::size_t... I>
auto fluid(std::index_sequence<I...>) -> State<Density, decltype((void) I, Vector{})..., Energy>;
}	 // namespace detail

/// One scalar followed by exactly ndim vector components.
template <typename Scalar, typename Vector>
using ScalarVectorState = decltype(detail::scalarVector<Scalar, Vector>(std::make_index_sequence<ndim>{}));

/// Density, exactly ndim momentum/velocity components, and energy/pressure.
template <typename Density, typename Vector, typename Energy>
using FluidState = decltype(detail::fluid<Density, Vector, Energy>(std::make_index_sequence<ndim>{}));

}	 // namespace octotigerII::units
